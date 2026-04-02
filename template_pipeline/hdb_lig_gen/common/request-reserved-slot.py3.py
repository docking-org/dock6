#!/usr/bin/env python
#
# request-reserved-slot.py
# Teague Sterling 2015
#
# Creates a simple, opt-in, process-based queue that can be 
# accessed by the file system. Intended to allow a process
# to request a specific core to limit to via taskset
# 
# Example:
# SLOT=./request-reserved-slot.py 
# taskset -c $SLOT some-problematic-command
#

import lockfile
import os
import psutil
import stat
import sys

WORLD_PERMISSIONS = stat.S_IRUSR | stat.S_IWUSR | \
                    stat.S_IRGRP | stat.S_IWGRP | \
                    stat.S_IROTH | stat.S_IWOTH
EXECUTE_PERMISSIONS = stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH

SYSTEM_LOCK_DIR = '/dev/shm/dock-process-reservations'
LOCK_FILE = os.path.join(SYSTEM_LOCK_DIR, 'lock')
RESERVATION_FILE = os.path.join(SYSTEM_LOCK_DIR, 'in-use')
LIMIT_FILE = os.path.join(SYSTEM_LOCK_DIR, 'max-reservations')


def get_processor_count():
    import multiprocessing
    return multiprocessing.cpu_count()


def setup():
    if not os.path.exists(SYSTEM_LOCK_DIR):
        os.makedirs(SYSTEM_LOCK_DIR)
        os.chmod(SYSTEM_LOCK_DIR, WORLD_PERMISSIONS | EXECUTE_PERMISSIONS)
    if not os.path.exists(RESERVATION_FILE):
        open(RESERVATION_FILE, 'w').close()
        os.chmod(RESERVATION_FILE, WORLD_PERMISSIONS)
    if not os.path.exists(LIMIT_FILE) or os.stat(LIMIT_FILE).st_size == 0:
        with open(LIMIT_FILE, 'w') as f:
            f.write("{0:d}".format(get_processor_count()))
        os.chmod(LIMIT_FILE, WORLD_PERMISSIONS)


def get_reservation_number(pid):
    with open(LIMIT_FILE, 'r') as f:
        max_reservations = int(f.read())
    slots = [None] * max_reservations
    with lockfile.FileLock(LOCK_FILE) as lock:
        with open(RESERVATION_FILE, 'r+') as f:
            for num, line in enumerate(f):
                line = line.strip()
                opid = int(line) if line else None
                slots[num] = opid
            if slots.count(None) == 0:
                for num, opid in enumerate(slots):
                    if opid is not None and not psutil.pid_exists(opid):
                        slots[num] = None
            try:
                reservation_number = slots.index(None)
                slots[reservation_number] = pid
                f.seek(0)
                f.truncate()
                f.write("\n".join('' if pid is None else str(pid) for pid in slots))
            except ValueError:
                reservation_number = None
    return reservation_number


def main(args, stdin=sys.stdin, stdout=sys.stdout, stderr=sys.stderr):
    setup()
    if len(args) > 0:
        reservation_pid = int(args[0])
    else:
        reservation_pid = os.getppid()
    reservation_number = get_reservation_number(reservation_pid)
    if reservation_number is not None:
        stdout.write("{0:d}\n".format(reservation_number))
        return 0
    else:
        stderr.write("No available reservations in {0}\n".format(SYSTEM_LOCK_DIR))
        return -1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

