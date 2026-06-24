import os
import subprocess
import sys
import time
from datetime import datetime

from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler

QUEUE_DIR = "analyzer_queue"
LOG_DIR = "logs"
MAX_PARALLEL_ANALYZERS = 3
ANALYZER_SCRIPT = "analyzers/post_ssmt.py"
#ANALYZER_SCRIPT = "analyzers/post_analysis.py"

os.makedirs(LOG_DIR, exist_ok=True)
os.makedirs(QUEUE_DIR, exist_ok=True)

running = {}

def log(msg):
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - {msg}")

def launch_analyzer(exp_id, exp_name, exp_type):
    """
    Lunch analyzer for experiment with subprocess in background
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
    log(f"Running analyzer for {exp_id}")
    log(f"Log for analyzer with {exp_id} at {str(log_path)}")
    log_file = open(log_path, "w", encoding="utf-8")
    os.environ["NO_COLOR"] = "1"  # disable color in subprocess's log
    try:
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
    except Exception as e:
        log(f"Failed to launch analyzer for {exp_id}: {e}")
        if log_file:
            log_file.close()
        return None, None


class ReadyFileHandler(FileSystemEventHandler):
    def on_created(self, event):
        """
        On Created events.
        Args:
            event: Event details.
        """
        log(f"[DEBUG] on_created triggered for {event.src_path}")
        if not event.is_directory and event.src_path.endswith('.ready'):
            if len(running) >= MAX_PARALLEL_ANALYZERS:
                log(f"Max parallel analyzers reached ({MAX_PARALLEL_ANALYZERS}), skipping {event.src_path}")
                return
            self.process_file(event.src_path)

    def process_file(self, filepath):
        base, _ = os.path.splitext(filepath)
        working_file = base + ".working"
        log(f"Detected .ready file: {filepath}")
        if os.path.exists(working_file):
            log(f"Working file already exists, skipping: {working_file}")
            return

        # Step 1: Rename to .working
        try:
            os.rename(filepath, working_file)
            log(f"Renamed to .working: {working_file}")
        except Exception as e:
            log(f"Failed to rename to .working: {e}")
            return

        # Step 2: Parse the flag content
        # Parse the flag content
        with open(working_file) as f:
            line = f.read().strip()
        parts = list(map(str.strip, line.split(',')))
        if len(parts) != 3 or not all(parts):
            log(f"Invalid format in {working_file}: {line}")
            return
        exp_id, exp_name, exp_type = parts

        # Step 3: Launch analyzer subprocess
        proc, log_file = launch_analyzer(exp_id, exp_name, exp_type)
        if not proc:
            return
        running[proc.pid] = (proc, log_file, working_file)


class Watcher:
    """
    Watches the bridge directory and communicates jobs to slurm.
    """

    def __init__(self, directory_to_watch: str, check_every: int = 5):
        """
        Creates our watcher.
        Args:
            directory_to_watch: Directory to sync from
            check_every: How often should the directory be synced
        """
        self.observer = Observer()
        self._directory_to_watch = directory_to_watch
        self._check_every = check_every

    def run(self):
        """
        Run the watcher.
        """
        event_handler = ReadyFileHandler()
        # Scans for pre-existing .ready files on startup
        for fname in os.listdir(self._directory_to_watch):
            fpath = os.path.join(self._directory_to_watch, fname)
            if fname.endswith('.ready') and os.path.isfile(fpath):
                log(f"[INIT SCAN] Found existing .ready file: {fpath}")
                event_handler.process_file(fpath)

        # Start watching the directory for new .ready files
        log(f"Watching {self._directory_to_watch} for .ready files")
        self.observer.schedule(event_handler, path=str(self._directory_to_watch), recursive=False)
        self.observer.start()
        try:
            while True:
                self.check_running_processes()
                time.sleep(self._check_every)
        except KeyboardInterrupt:
            log("KeyboardInterrupt received. Stopping observer.")
            self.observer.stop()
        except Exception as e:
            log(f"Unexpected error: {e}")
            self.observer.stop()
        self.observer.join()

    def check_running_processes(self):
        """
        Check running processes to see if any analyzers have completed.
        If so, rename the .working file to .done or .fail. Then delete the process id from the running dict.
        Returns:

        """
        if running:
            log(f"Running analyzers: {len(running)}")
        finished = []
        for pid, (proc, log_file, flag_path) in running.items():
            if proc.poll() is not None:
                exit_code = proc.returncode
                log_file.close()
                finished.append(pid)
                if exit_code == 0:
                    done_path = flag_path.replace(".working", ".done")
                    os.rename(flag_path, done_path)
                    log(f"Analyzer completed: {flag_path} → {done_path}")
                else:
                    fail_path = flag_path.replace(".working", ".fail")
                    os.rename(flag_path, fail_path)
                    log(f"Analyzer failed: {flag_path} → {fail_path}")

        for pid in finished:
            del running[pid]


def main():
    w = Watcher(directory_to_watch=QUEUE_DIR)
    w.run()

if __name__ == "__main__":
    log("Analyzer manager running...")
    main()
