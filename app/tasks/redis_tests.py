import time

import logging
logger = logging.getLogger()
def example(seconds):
    seconds=int(seconds)
    for i in range(seconds):
        logger.debug(i)
        time.sleep(1)
    logging.debug("Done")
# example(10)
