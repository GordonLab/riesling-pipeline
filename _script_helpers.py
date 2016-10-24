"""
Shared helper functions.

License: MIT
Author: Nick Semenkovich <semenko@alum.mit.edu>
"""

from __future__ import print_function
import glob
import os
import shlex
import time
import tempfile
import yaml
from distutils import spawn
from subprocess import Popen, PIPE


if __name__ == '__main__':
    print("Do not run this as a standalone script.")
    exit()

###################
## Early sanity checks.
###################

THISPATH = os.path.dirname(os.path.realpath(__file__))

## TODO: Given our new dist, consider dropping some of these or modifying this.
for cmd in ['samtools', 'grep', 'bedtools', 'bowtie2']:
    if spawn.find_executable(cmd) is None:
        raise OSError("Software missing, unable to find: %s" % (cmd))


def get_config():
    """
    Load the yaml config based on our true root path.

    :return: YAML config object
    """
    with open(THISPATH + '/.config.yaml') as yamlfile:
        config = yaml.load(yamlfile)

    # Append our true path to each binary
    binaries_with_paths = {}
    for binary_name in config['binaries']:
        binaries_with_paths[binary_name] = THISPATH + '/' + config['binaries'][binary_name]

    config['binaries'] = binaries_with_paths

    return config

def setup_output_path(path_or_file):
    """
    Make sure our output directory is writeable. Create it if necessary.
    """
    if os.path.isfile(path_or_file):
        raise ValueError("Output path appears to be a file. Please specify a directory.")

    output_path = path_or_file

    output_path = os.path.normpath(os.path.normcase(output_path))

    try:
        os.mkdir(output_path)
    except OSError:
        if not os.access(output_path, os.W_OK):
            raise OSError("Output path couldn't be created or isn't writeable: %s" % (output_path))

    return output_path


def validate_input_files(input_path_or_file, mask='.bam'):
    """
    Given an input arg (either a specific file, or a path), return it as a list of files.
    Also check that files are readable.
    """
    # TODO: Make this handle a list of raw files, too. (e.g. for 3-OPTIONAL...)
    if os.path.isfile(input_path_or_file):
        # We got a single file!
        if not input_path_or_file.endswith(mask):
            raise ValueError("Expected a %s input (or a directory with %s). You gave: %s"
                             % (mask, mask, input_path_or_file))
        file_list = [input_path_or_file]
    else:
        # It's not a file. Must be a directory.
        if not os.path.isdir(input_path_or_file):
            raise ValueError("Input not found (or not a file/folder): %s" % ((input_path_or_file)))
        file_list = glob.glob((input_path_or_file) + "/*" + mask)

    if len(file_list) == 0:
        raise ValueError("Input was empty!")

    for filename in file_list:
        if not os.access(filename, os.R_OK):
            raise OSError("Cannot read file: %s" % (filename))

    return file_list


class ShellJobRunner():
    """
    Run shell jobs and make sure they complete.

    This is dangerous to run on untrusted inputs!
    """
    def __init__(self, logger, delay_seconds=False):
        self.logger = logger
        self.delay_seconds = delay_seconds
        self.process_list = []
        if delay_seconds is False:
            self.logger.info('Created a NON-parallel job runner.')
        else:
            self.logger.info('Created a parallel job runner with %i second delay between jobs.' % (delay_seconds))

    def run(self, command):
        """
        Run a given command. May be blocking (default) or non-blocking if delay_seconds is set.
        """

        self.logger.debug('Running: %s' % (command))

        # TODO: Clean this up? We shouldn't be spawning sh to spawn bash to set pipefail ...
        process = Popen('nice bash -c "set -o pipefail; (%s)"' % command, shell=True)
        self.logger.debug('Spawned PID: %i' % (process.pid))
        self.process_list.append(process)

        if self.delay_seconds is False:
            self.logger.info('* Parallelism disabled. Waiting for job to complete.')
            runtime_process_status = process.wait()
        else:
            self.logger.info('* Waiting %i seconds to spawn next job.' % (self.delay_seconds))
            time.sleep(self.delay_seconds)
            runtime_process_status = process.poll()

            if runtime_process_status is None:
                # Not done yet, that's cool!
                # Since delay_seconds was set, we'll return now. The user better call finish() later!
                pass
            elif runtime_process_status == 0:
                # We're done already? That was suspiciously fast (or delay_seconds is too high).
                self.logger.warn('This task finished in less than %d seconds.' % (self.delay_seconds))
                self.logger.warn('This is OK if your input files are small, otherwise, this is suspicious.')

        if runtime_process_status > 0:
            self.logger.critical('The last command failed!')
            self.logger.critical('Fault occurred in: %s' % (command))
            raise ValueError('Process failed with exit code: %i' % (runtime_process_status))


    def finish(self):
        """
        Close out / block for processes.
        """
        self.logger.info('Waiting for all %i processes to complete...' % (len(self.process_list)))

        # TODO: Consider more granular failure info here?
        exit_codes = [p.wait() for p in self.process_list]
        if sum(exit_codes) != 0:
            self.logger.critical('A process died! Cannot continue.')
            raise ValueError("One of the processes failed! Are you out of RAM (or hitting a system limit?)")

        self.logger.info('All processes done! Yay!')


class IntelligentRunner():
    """
    Run the input command string (echo | grep | cut ...) via subprocess, and
    catch / discard known false-positive/annoying errors.
    """
    known_ignorable_stderr = {'[samopen] SAM header is present:'}
    stderr_fp = tempfile.SpooledTemporaryFile()

    def __init__(self, command_string):
        self.command_string = command_string
        self.command_list = ','.join(shlex.split(command_string)).split(",|,")

    def _check_for_errors(self):
        self.stderr_fp.seek(0)
        for stdout_line in self.stderr_fp.readlines():
            if stdout_line.strip() not in self.known_ignorable_stderr:
                if len(filter(stdout_line.strip().startswith, self.known_ignorable_stderr)) == 0:
                    print("Was running: %s" % (self.command_string))
                    raise ValueError("Unignorable STDERR output: %s" % (stdout_line.strip()))
            print(stdout_line.strip())


    def run(self):
        print(self.command_list)
        last_process = Popen(self.command_list[0].split(','), stdout=PIPE, stderr=self.stderr_fp)
        for command in self.command_list:
            print(command)
            last_process = Popen(command.split(','), stdin=last_process.stdout, stdout=PIPE, stderr=self.stderr_fp)

        # Grab the output
        output = last_process.communicate()[0]

        self._check_for_errors()

        return output
