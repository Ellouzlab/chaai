import os
import logging

from chaai.cli.parser import argparser
from chaai.cli.logsetup import init_logging


_COMMAND_HANDLERS = {
    'build_db': ('chaai.cli.commands.build_db', 'BuildDbHandler'),
    'dist':    ('chaai.cli.commands.dist',     'DistHandler'),
    'search':  ('chaai.cli.commands.search',   'SearchHandler'),
}


def _dispatch(command: str, arguments):
    from importlib import import_module

    module_path, class_name = _COMMAND_HANDLERS[command]
    module = import_module(module_path)
    handler_cls = getattr(module, class_name)
    handler_cls(arguments)


def main():
    arguments = argparser()

    log_path = os.path.join("chaai_logs")
    os.makedirs(log_path, exist_ok=True)

    command = arguments.command
    init_logging(os.path.join(log_path, f"{command}.log"))
    logging.info(f"chaai {command}")
    _dispatch(command, arguments)


if __name__ == "__main__":
    main()
