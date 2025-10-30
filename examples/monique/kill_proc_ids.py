import os
import signal

PID_FILE = "pid_list.txt"

if not os.path.exists(PID_FILE):
    print(f"{PID_FILE} not found.")
    exit(1)

with open(PID_FILE) as f:
    pids = [int(line.strip()) for line in f if line.strip().isdigit()]

for pid in pids:
    try:
        # On Windows, SIGTERM translates to a termination signal.
        os.kill(pid, signal.SIGTERM)
        print(f"Killed PID {pid}")
    except ProcessLookupError:
        print(f"PID {pid} not found (already exited).")
    except PermissionError:
        print(f"No permission to kill PID {pid}.")

os.remove(PID_FILE)
print(f"Removed {PID_FILE}")