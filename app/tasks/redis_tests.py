import time
import logging


logging.basicConfig(level=logging.DEBUG)

def example(seconds):
    logging.info("************ Starting task *****************")

    for i in range(seconds):
        logging.debug(i)
        time.sleep(1)
    logging.info("************ Task complete *****************")
