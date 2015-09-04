import os
import re
from subprocess import Popen, STDOUT, PIPE
import sys

'''
This script is a wrapper around scons to facilitate the build
process of the lapack_functions. It tries to build, and if there
is a failure due to sconsflags set in the SCONSFLAGS environment 
variable the offending flag is removed, and the build re-tried.
This will happen up to a maximum of 100 tries, after which it
will print a failure message and exit. If the build fails due to
any other reason than a flag error, the script will print the
output and exit.
'''

def printError(output):
    print("There was an error in the build, please see")
    print(output)
    sys.exit(1)

# Get a copy of the current environment
thisEnv = os.environ.copy()

# Watch the exitcode in a loop, once the build succeeds, exit
exitcode = 2
# While loops can be dangerous, as they can loop forever. As such
# create a variable that if the loop exceeds the limit it will break
failsafe = 100
counter = 0
# Let's keep track of anything which is not used
notUsed = []

while exitcode > 0:
    counter += 1
    process = Popen(['scons', '-Q'], stderr=STDOUT,
                    stdout=PIPE, env=thisEnv)
    # Wait for the build to finish and get the return code
    rcode = process.wait()
    # Read the output (as well as stderr)
    output = process.stdout.read()
    # If the script succeeded break the loop and print the output
    if rcode == 0:
        exitcode = 0
        # If we have thrown out some flags, let the user know
        if len(notUsed) > 0:
            print('The following flags caused the build to fail' +
                  'and were not used')
            print(notUsed)
        print(output)
    else:
        # If the counter exceeds the limit, fail out immediately
        if counter > failsafe:
            print("The build could not scceed")
            printError(output)
        # There was and issue with the command, process the output
        exitcode = rcode
        # Check if SCONSFLAGS exists and is no zero length, if either
        # of these are not true, then the build failed for some reason
        # not related to flags and we will handle that below.
        if "SCONSFLAGS" in thisEnv and\
                len(thisEnv['SCONSFLAGS']) > 0:
            # Get the offending sconsflag key
            problemKey = re.findall('.*option: (.*)\n*', output)[0]
            # Check if we were unable to extract the key, or if the key
            # was not in the environment variable, if so print error
            # and exit
            if len(problemKey) == 0 or \
               problemKey not in thisEnv['SCONSFLAGS']:
                printError(output)
            else:
		thisEnv['SCONSFLAGS'] = thisEnv['SCONSFLAGS']\
                                        .replace(problemKey,'')
                notUsed.append(problemKey)
        else:
            # The error does not seem to be related to flags
            # print the output and exit
            printError(output)
