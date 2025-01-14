from js_wrappers.filesystem import DirectoryHandle
from js_wrappers.async_js import wrap_promise, sleep_a
from js import setWorkingDirectory, getWorkingDirectory

async def set_working_directory_a():
    handle = DirectoryHandle()
    await handle.open_a()
    if handle.is_open():
        setWorkingDirectory(handle._handle)

async def get_working_directory_a():
    handle = await wrap_promise(getWorkingDirectory())
    if not handle:
        return None
    result = DirectoryHandle()
    result._handle = handle
    return result

async def prompt_get_working_directory_a():
    print("Will prompt to set working directory momentarily.")
    await sleep_a(1000)
    print("Choice will persist until you clear brower memory or reset working directory.")
    await sleep_a(1000)
    await set_working_directory_a()
