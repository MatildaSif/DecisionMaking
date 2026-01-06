#!/usr/bin/env bash

# Parameter Recovery Parallel Execution Script
# Launches multiple tmux panes running parameter recovery jobs

SESSION_NAME="param_recovery"
NUM_JOBS=20  # Number of parallel jobs to run
SCRIPT_NAME="param_recovery_batch.R"

# Check if the session already exists
if tmux has-session -t $SESSION_NAME 2>/dev/null; then
    echo "Session $SESSION_NAME already exists. Attaching to it."
    tmux attach-session -t $SESSION_NAME
    exit 0
fi

# Create new session with first pane
echo "Creating tmux session: $SESSION_NAME"
tmux new-session -d -s $SESSION_NAME

# Launch first job in the initial pane
tmux send-keys -t $SESSION_NAME:0 "Rscript $SCRIPT_NAME 1000" C-m

# Create additional panes and launch jobs
for i in $(seq 1 $((NUM_JOBS - 1))); do
    # Calculate seed offset for this job
    SEED_OFFSET=$((1000 + i))
    
    # Split window vertically
    tmux split-window -t $SESSION_NAME -h
    
    # Rebalance panes (optional, helps with layout)
    if [ $((i % 5)) -eq 0 ]; then
        tmux select-layout -t $SESSION_NAME tiled
    fi
    
    # Send command to the new pane
    tmux send-keys -t $SESSION_NAME "Rscript $SCRIPT_NAME $SEED_OFFSET" C-m
    
    echo "Launched job $i with seed offset $SEED_OFFSET"
done

# Set optimal layout
tmux select-layout -t $SESSION_NAME tiled

# Attach to the session
echo "Attaching to tmux session $SESSION_NAME"
echo "Press Ctrl+B then D to detach, or Ctrl+B then X to kill a pane"
tmux attach-session -t $SESSION_NAME