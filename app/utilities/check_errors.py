import sys
import os
import logging
logger = logging.getLogger('routes')


def count_lines(fname):
    """ Count number of lines in file.

    :param fname: Name (path) of the file to be counted.
    :type  fname: str (path)

    """

    try:
        with open(fname, 'r') as fh:
            number_of_lines = len(fh.readlines())

    except IOError:
        logging.warning("%s is empty.", fname)

        return 0

    except:
        raise

    return number_of_lines


def rayt_repin_counts(outdir,
                      reference_strain,
                      ):
    """ Parse output files from a given job and infer how the job has ended.

    :param outdir: The rarefan run output directory.
    :type  run_id: str (path)

    :param reference_strain: The reference strain
    :type  reference_strain: str

    :return: A dictionary of counts for number of rayts, number of overrepresented reps and number of repins.
    :rtype: dict

    """

    number_of_rayts = count_lines(os.path.join(outdir, "repin_rayt_association.txt.fas")) // 2
    number_of_overreps = count_lines(os.path.join(outdir, reference_strain + ".overrep"))

    results = {"number_of_rayts": number_of_rayts,
               "number_of_overreps": number_of_overreps,
               "number_of_repins": {}
               }

    for i in range(0, 5):
        line_count = count_lines(
            os.path.join(
                outdir,
                f"reference_strain_{i}",
                f"reference_strain_{i}_largestCluster.nodes"
            )
        )

        results["number_of_repins"][i] = line_count

    return results
