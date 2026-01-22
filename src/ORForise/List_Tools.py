from importlib import import_module
import  os

try:
    from .utils import *
except (ImportError, ModuleNotFoundError):
    from utils import *





def main():
    print('ORForise ' + ORForise_Version + ': List Tools Run Parameters')

    tools = set()
    base_dirs = [
        os.path.join(os.path.dirname(__file__), 'Tools'),
        os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Tools')),
    ]

    for base in base_dirs:
        if not os.path.isdir(base):
            continue
        try:
            for entry in os.listdir(base):
                entry_path = os.path.join(base, entry)
                if os.path.isdir(entry_path) and not entry.startswith('.') and entry != '__pycache__':
                    tools.add(entry)
        except OSError:
            continue

    if not tools:
        print('No tools found in the searched directories.')
        return

    print('Available tools:')
    for tool_name in sorted(tools):
        print(' -', tool_name)
        try:
            tool_ = import_module('Tools.' + tool_name + '.' + tool_name)
            print('    Imported from Tools.' + tool_name)
        except ModuleNotFoundError:
            try:
                tool_ = import_module('ORForise.Tools.' + tool_name + '.' + tool_name)
                print('    Imported from ORForise.Tools.' + tool_name)
            except ModuleNotFoundError:
                print('    Tool not importable')



if __name__ == "__main__":
    try:
        try:
            main()
        except Exception:
            print('Unhandled exception in main')
    finally:
        print(CLOSING)