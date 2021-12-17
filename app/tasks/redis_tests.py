# usage: $> rq enqueue -p app.tasks.redis_tests example 10
import time

import logging
logger = logging.getLogger(__name__)

def example(seconds):
    seconds=int(seconds)
    for i in range(seconds):
        logger.debug(i)
        time.sleep(1)
    logger.debug("Done")
