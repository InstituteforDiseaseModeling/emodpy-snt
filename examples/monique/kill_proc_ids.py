import os
import signal
import platform

PID_FILE = "pid_list.txt"

if not os.path.exists(PID_FILE):
    print(f"{PID_FILE} not found.")
    exit(1)

with open(PID_FILE) as f:
    pids = [int(line.strip()) for line in f if line.strip().isdigit()]

for pid in pids:
    try:
        if platform.system() == "Windows":
            os.system(f"taskkill /PID {pid} /F >nul 2>&1")
        else:
            os.kill(int(pid), signal.SIGTERM)
        print(f"Killed PID {pid}")
    except Exception as e:
        print(f"Failed to kill PID {pid}: {e}")

os.remove(PID_FILE)
print(f"Removed {PID_FILE}")