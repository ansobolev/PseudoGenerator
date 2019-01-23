
import logging

def get_logger(name, element):
    # create logger
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    # create console handler and set level to debug
    fh = logging.FileHandler(element + '/log.dat')
    fh.setLevel(logging.DEBUG)

    # create formatter
    formatter = logging.Formatter('%(asctime)s: %(name)s - %(levelname)s - %(message)s')

    # add formatter to fh
    fh.setFormatter(formatter)

    # add fh to logger
    logger.addHandler(fh)
    return logger

def interlog(logger):
    logger.info('-----' * 10)
