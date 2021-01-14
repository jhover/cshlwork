import os


def get_taskid():
    tid = 1
    try:
        tid = int(os.environ['SGE_TASK_ID'])
    except KeyError:
        pass
    return tid


