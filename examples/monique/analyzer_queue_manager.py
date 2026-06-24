# Another way to run analyzer separately from run_simulation.py with centralized queue manager, so you can:
# Control how many analyzers run at once
# Prevent memory/CPU exhaustion
# Cleanly separate experiment submission from analysis
import os
import signal
import time
import subprocess
import sys
from datetime import datetime

QUEUE_DIR = "analyzer_queue"
LOG_DIR = "logs"
MAX_PARALLEL_ANALYZERS = 2
#ANALYZER_SCRIPT = "analyzers/post_ssmt.py"
ANALYZER_SCRIPT = "analyzers/post_analysis.py"

running = {}

def log(msg):
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - {msg}")

def launch_analyzer(exp_id, exp_name, exp_type):
    """
    lunch analyzer for experiment with subprocess in background
    Args:
        exp_id: experiment id, str
        exp_name: experiment name, str
        exp_type: experiment type, str

    Returns:

    """
    if "ssmt" in ANALYZER_SCRIPT:
        log_path = os.path.join(LOG_DIR, f"ssmt_{exp_name}_{exp_id}.log")
    else:
        log_path = os.path.join(LOG_DIR, f"local_{exp_name}_{exp_id}.log")
    print(f"Running analyzer for {exp_id}")
    print(f"Log for analyzer with {exp_id} at {str(log_path)}")
    log_file = open(log_path, "w", encoding="utf-8")
    os.environ["NO_COLOR"] = "1"  # disable color in subprocess's log
    proc = subprocess.Popen(
        [
            sys.executable, ANALYZER_SCRIPT,
            "--exp-id", exp_id,
            "--type", exp_type,
            "--name", exp_name
        ],
        stdout=log_file,
        stderr=subprocess.STDOUT,
        text=True
    )
    return proc, log_file

def graceful_exit(signum, frame):
    print("Stopping analyzer manager...")
    sys.exit(0)


def main():
    """
    Consistently run analyzers until interrupted.
    Only pick files from the queue directory with .ready extension.
    Returns:

    """
    signal.signal(signal.SIGINT, graceful_exit)
    os.makedirs(QUEUE_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)

    # run consistently until interrupted
    while True:
        # Check for new .ready files
        ready_files = [f for f in os.listdir(QUEUE_DIR) if f.endswith(".ready")]
        for fname in ready_files:
            if len(running) >= MAX_PARALLEL_ANALYZERS:
                break
            if not fname.endswith(".ready"):
                continue

            ready_path = os.path.join(QUEUE_DIR, fname)
            working_path = ready_path.replace(".ready", ".working")
            print(f"Replace: {ready_path} → {working_path}")

            try:
                os.rename(ready_path, working_path)
            except FileNotFoundError:
                continue  # Another process already picked it up

            # Parse the flag content
            with open(working_path) as f:
                line = f.read().strip()
                exp_id, exp_name, exp_type = line.split(',')

            proc, log_file = launch_analyzer(exp_id, exp_name, exp_type)
            running[proc.pid] = (proc, log_file, working_path)

        # Check running processes to see if any analyzers have completed
        finished = []
        for pid, (proc, log_file, flag_path) in running.items():
            # If the process has finished (not None)
            if proc.poll() is not None:
                log_file.close()
                done_path = flag_path.replace(".working", ".done")
                os.rename(flag_path, done_path)
                print(f"Analyzer completed: {flag_path} → {done_path}")
                finished.append(pid)

        for pid in finished:
            del running[pid]

        time.sleep(5)

if __name__ == "__main__":
    print("Analyzer manager running...")
    main()
