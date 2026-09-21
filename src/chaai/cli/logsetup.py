import os
import sys
import logging


def get_unique_log_filename(base_log_filename: str) -> str:
    i = 1
    base_name, extension = os.path.splitext(base_log_filename)
    while os.path.exists(base_log_filename):
        base_log_filename = f"{base_name}_{i}{extension}"
        i += 1
    return base_log_filename


class TeeHandler(logging.Handler):

    def __init__(self, filename: str, mode: str = 'a'):
        super().__init__()
        self.file = open(filename, mode)
        self.stream_handler = logging.StreamHandler(sys.stdout)

    def emit(self, record):
        log_entry = self.format(record)
        self.file.write(f"{log_entry}\n")
        self.file.flush()
        self.stream_handler.emit(record)

    def close(self):
        self.file.close()
        super().close()


def init_logging(log_filename: str) -> None:
    while os.path.exists(log_filename):
        log_filename = get_unique_log_filename(log_filename)
    print(f"Logging to {log_filename}")

    logger = logging.getLogger()

    for h in list(logger.handlers):
        logger.removeHandler(h)

    logger.setLevel(logging.INFO)

    tee_handler = TeeHandler(log_filename)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    tee_handler.setFormatter(formatter)

    logger.addHandler(tee_handler)
