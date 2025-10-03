import os
import json
import subprocess
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Run a pipeline based on a configuration file.")
parser.add_argument("--config", required=True, help="Path to the configuration file (e.g., config_step1.json)")
args = parser.parse_args()

# Load configuration file
config_path = args.config
with open(config_path, 'r') as config_file:
    config = json.load(config_file)

# Extract the base paths
script_base_path = config["script_base_path"]
data_base_path = config["data_base_path"]

# Define the pipeline steps dynamically from the config
pipeline = config["pipeline"]

# Run the pipeline
for step in pipeline:
    script = os.path.join(script_base_path, step["script"])  # Prepend script_base_path
    step_config = config["steps"][step["config_key"]]

    # Build the command dynamically
    cmd = ["python", script]
    if isinstance(step_config["input"], list):
        cmd.extend([os.path.join(data_base_path, path) for path in step_config["input"]])
    else:
        cmd.append(os.path.join(data_base_path, step_config["input"]))
    if isinstance(step_config["output"], list):
        cmd.extend([os.path.join(data_base_path, path) for path in step_config["output"]])
    else:
        cmd.append(os.path.join(data_base_path, step_config["output"]))

    # Print and execute the command
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Check for errors
    if result.returncode != 0:
        print(f"Error in {script}: {result.stderr}")
        break
    else:
        print(f"Completed: {script}")
        print(result.stdout)

print("Pipeline completed.")