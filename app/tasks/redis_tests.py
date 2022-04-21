# usage: $> rq enqueue -p app.tasks.redis_tests example 10
import time

from app import app
app.app_context().push()

logger = app.logger

def example(seconds):
    seconds=int(seconds)
    for i in range(seconds):
        logger.debug(i)
        time.sleep(1)
    logger.debug("Done")
