#!/usr/bin/env python3
"""List available ORForise tools (very simple, no CLI).

For each entry under a discovered Tools directory this script will try to import
Tools.<tool>.<tool> and then ORForise.Tools.<tool>.<tool>. It prints a tab-separated
list: <tool_name>\tOK/NOT_AVAILABLE\t<note>
"""

from importlib import import_module
import os
import sys


def find_tools_dir():
    # Try several likely locations and return the first existing directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.abspath(os.path.join(script_dir, '..', '..'))
    candidates = [
        os.path.join(os.getcwd(), 'Tools'),
        os.path.join(os.getcwd(), 'tools'),
        os.path.join(repo_root, 'Tools'),
        os.path.join(repo_root, 'tools'),
        os.path.join(repo_root, 'ORForise', 'Tools'),
        os.path.join(repo_root, 'galaxy', 'tools'),
    ]
    for c in candidates:
        if os.path.isdir(c):
            return c
    return None


def main():
    tools_dir = find_tools_dir()
    if not tools_dir:
        print('Tools directory not found')
        return 1

    entries = sorted([e for e in os.listdir(tools_dir) if not e.startswith('.')])
    for name in entries:
        # Attempt to import as a tool module using the project's two-location pattern
        try:
            import_module('Tools.' + name + '.' + name, package='my_current_pkg')
            print(f"{name}\tOK\timported from Tools.{name}")
            continue
        except ModuleNotFoundError:
            pass

        try:
            import_module('ORForise.Tools.' + name + '.' + name, package='my_current_pkg')
            print(f"{name}\tOK\timported from ORForise.Tools.{name}")
            continue
        except ModuleNotFoundError:
            pass

        # Not importable as a tool module
        print(f"{name}\tNOT_AVAILABLE")

    return 0


if __name__ == '__main__':
    sys.exit(main())
