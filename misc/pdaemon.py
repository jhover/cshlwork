#!/usr/bin/env python
import daemon
import os
import time

LOGFILE=os.path.expanduser("~/pdaemon.log")

def my_daemon_program():
    cycles = 0
    with open(LOGFILE, 'a', encoding = 'utf-8') as f:
        while True:
            f.write(f"I'm a daemon cycle={cycles}\n")
            f.flush()
            cycles += 1
            time.sleep(15)

with daemon.DaemonContext():
    my_daemon_program()    
    