// import Worker from './pyodide.worker.js';
import { v4 as uuid4 } from "uuid";
import { sleep } from "./utils";
import { EventEmitter } from "eventemitter3";

export class PythonExecutor {
    constructor(){
        this.executions = {};
        this.completers = {};
        this.pyodide_worker = new Worker("pyodide_worker.bundle.js");
        this.pyodide_worker.addEventListener("message", this._handleMessage.bind(this));
        // The pyodide worker needs to be able to send messages to the service worker, so we make a channel
        // and send one end to the service worker and the other to the pyodide worker.
        let {port1, port2} = new MessageChannel();
        this.pyodide_worker.postMessage({
            cmd : "service_worker_channel", 
            port : port1,
        }, [port1]);
        navigator.serviceWorker.controller.postMessage({
            cmd : "pyodide_worker_channel", 
            port : port2,
        }, [port2]);

        let _readyPromise = new Promise((resolve, reject) => this._readyPromise = {resolve, reject});
        this._readyPromise.promise = _readyPromise;
    }

    _handleMessage(event){
        let message = event.data;
        let message_cmd = message.cmd;
        let subhandler_name = ({"execute" : "_handleExecutionMessage", "complete" : "_handleCompletionMessage", "ready" : "_handleReadyMessage"})[message_cmd];
        if(!subhandler_name){
            throw new Error(`Unknown command "${message_cmd}"`);
        }
        this[subhandler_name](message);
    }
    
    _handleReadyMessage(message){
        if(message.exception){
            this._readyPromise.reject(message.exception);
            return;
        }
        console.log("Pyodide is ready!");
        this._readyPromise.resolve();
    }

    _handleExecutionMessage(message){
        // execution messages get emitted on the execution object.
        const { uuid, subcmd, last_response } = message;
        const execution = this.executions[uuid];
        if(!execution){
            throw new Error(`Invalid execution uuid "${uuid}"`);
        }
        // Check if there is a handler for the given command, otherwise fail.
        // All messages are meant to be handled.
        if(execution.listenerCount(subcmd) === 0) {
            throw new Error(`Unexpected command "${subcmd}"`);
        }
        execution.emit(subcmd, message);
        if(last_response){
            execution._close();
            delete this.executions[uuid];
        }
    }

    _handleCompletionMessage(message){
        const { uuid, subcmd } = message;
        const completer = this.completers[uuid];
        if(!completer){
            throw new Error(`Invalid completer uuid "${uuid}"`);
        }
        if(completer.listenerCount(subcmd) === 0) {
            throw new Error(`Unexpected command "${subcmd}"`);
        }
        completer.emit(subcmd, message);
    }

    _postMessage(cmd, uuid, msg){
        Object.assign(msg, {cmd, uuid});
        this.pyodide_worker.postMessage(msg);
    }

    async ready(){
        return await this._readyPromise.promise;
    }


    execute(code){
        const interrupt_buffer = new Int32Array(new SharedArrayBuffer(4));
        const uuid = uuid4();
        const execution = new Execution(interrupt_buffer);
        this.executions[uuid] = execution;
        this._postMessage("execute", uuid, {code, interrupt_buffer});
        return execution;
    }

    new_completer(){
        const uuid = uuid4();
        const completer = new Completer(this, uuid);
        this.completers[uuid] = completer;
        this._postMessage("complete", uuid, {subcmd : "new_completer"});
        return completer;
    }

}

export class Execution extends EventEmitter {
    /* An execution object. This is for attaching handlers / giving out promises for various lifecycle events of the execution.
       The execution object is created and scheduled by PythonExecutor.execute. Other files do not construct these directly.
       The Executor also dispatches messages from the pyodide worker to the appropriate execution.
       See the python file "execution.py" for when python generates the messages this is responding to.
    */
    constructor(interrupt_buffer){
        super();
        this.interrupt_buffer = interrupt_buffer;
        // Using "once" here helps us throw a useful error if some logic error causes the pyodide worker to send 
        // the same event twice.
        this._validate_syntax = new Promise((resolve, reject) => {
            this.once("validate_syntax", resolve);
        });
        this._result = new Promise((resolve, reject) => {
            this.once("result", (message) => resolve(message.result));
            this.once("exception", (message) => reject(message));
            this.once("keyboard_interrupt", (message) => reject(message));
        });
    }
    
    async validate_syntax(){
        return await this._validate_syntax;
    }

    async result(){
        return await this._result;
    }

    setInterrupt(i){
        this.interrupt_buffer[0] = i;
        // Atomics.notify(this.interrupt_buffer, 0);
    }

    keyboardInterrupt(){
        this.setInterrupt(2); // SIGINT
    }

    onStdout(handler, context){
        this.on("stdout", function(message) { 
            handler.call(this, message.data); 
        }, context);
    }

    ignoreStdout(){
        this.on("stdout", () => undefined);
    }

    onStderr(handler, context){
        this.on("stderr", function(message) { 
            handler.call(this, message.data); 
        }, context);
    }

    ignoreStderr(){
        this.on("stderr", () => undefined);
    }

    _close(){
        
    }    
}

export class Completer extends EventEmitter {
    constructor(executor, uuid){
        super();
        this.executor = executor;
        this.uuid = uuid;
        this.responses = {};
        for(let cmd of ["completions", "completion_detail"]){
            this._attachResponseHandler(cmd);
        }
    }

    _attachResponseHandler(cmd){
        this.on(cmd, (msg) => {
            let promise_obj = this.responses[msg.subuuid];
            if(!promise_obj){
                throw Error(`Unknown subuuid ${subuuid}`);
            }
            if(cmd !== promise_obj.cmd) {
                throw new Error(`Wrong command for response subuuid ${subuuid}. Was expecting command to be "${cmd}" but got "${promise_obj.cmd}"`);
            }
            promise_obj.resolve(msg);
        });
    }

    _postMessage(subcmd, msg){
        Object.assign(msg, {subcmd});
        this.executor._postMessage("complete", this.uuid, msg);
    }

    async getCompletions(code, position){
        let subuuid = uuid4();
        let response_promise = new Promise((resolve, reject) => 
            this.responses[subuuid] = {resolve, reject, cmd : "completions"}
        );
        let {lineNumber, column} = position;
        this._postMessage("completions", { subuuid, code, lineNumber, column });
        let response = await response_promise;
        return [response.state_id, response.completions];
    }

    async getCompletionInfo(state_id, idx){
        let subuuid = uuid4();
        let response_promise = new Promise((resolve, reject) => 
            this.responses[subuuid] = {resolve, reject, cmd : "completion_detail"}
        );
        this._postMessage("completion_detail", { subuuid, idx, state_id });
        let {docstring, signature} = await response_promise;
        return {docstring, signature};
    }

    close(){
        delete this.executor.completers[this.uuid];
        delete this.executor;
    }

}