#!/usr/bin/env bash

##### HUGE KUDOS TO Balázs FOR SETTING THIS UP!!!! #####

#install packages
source setup.sh

SESSION_NAME="TMUX3"


# Check if the session already exists
if tmux has-session -t $SESSION_NAME 2>/dev/null; then
    echo "Session $SESSION_NAME already exists. Attaching to it."
    tmux attach-session -t $SESSION_NAME
else
    # Create a new session and panes 
    # (there must be a more automated way with a loop)
    tmux new-session -d -s $SESSION_NAME
    #tmux split-window -h
    #tmux split-window -h
    #tmux select-pane -t 0
    #tmux split-window -h
    #tmux select-pane -t 0
    #tmux split-window -v
    #tmux select-pane -t 2
    #tmux split-window -v
    #tmux select-pane -t 4
    #tmux split-window -v
    #tmux select-pane -t 6
    #tmux split-window -v
    
    for i in {0..20}
    do
      tmux split-window -v -p 5
      tmux select-pane -t $i
    done
    
    #loop for seeds and send them to the script
    for i in {0..20}
    do
      #define seed
      SEED=$((i + 1001))
      
      # build commands
      #tmux_command+="\; split-window $command"
      # Send a command to the pane, make R read SEED
      #tmux split-window
      #tmux select-pane -t $i
      tmux send-keys -t $i 'Rscript Tmux_ORL_rec.R ' $SEED C-m
    
    done
    #tmux "$SESSION_NAME attach"
 # Attach to the created session
    tmux attach-session -t $SESSION_NAME
fi