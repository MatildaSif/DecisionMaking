#!/usr/bin/env bash

##### HUGE KUDOS TO Balázs FOR SETTING THIS UP!!!! #####

#install packages
source setup.sh

SESSION_NAME="TMUX"


# Check if the session already exists
if tmux has-session -t $SESSION_NAME 2>/dev/null; then
    echo "Session $SESSION_NAME already exists. Attaching to it."
    tmux attach-session -t $SESSION_NAME
else
    # Create a new session and panes 
    # (there must be a more automated way with a loop)
    tmux new-session -d -s $SESSION_NAME
    
    # setting up 20 panes (4 columns by 5 rows)
    #tmux split-window -h -p 75
    #tmux split-window -h -p 66.67
    #tmux split-window -h -p 50
    #for i in {0..3}
    #do
    #  tmux select-pane -t $i
    #  tmux split-window -v -p 80
    #  tmux split-window -v -p 75
    #  tmux split-window -h -p 66.67
    #  tmux split-window -h -p 50
    #done
    
    tmux split-window -h -l 60
    tmux split-window -h -l 40
    tmux split-window -h -l 20
    for i in {0..3}
    do
      tmux select-pane -t $i
      tmux split-window -v -l 40
      tmux split-window -v -l 30
      tmux split-window -h -l 20
      tmux split-window -h -l 10
    done
    
    for i in {0..20}
    do
      
    done
    
    #loop for seeds and send them to the script
    for i in {0..20}
    do
      #define seed
      SEED=$((i + 1001))
      
      tmux select-pane -t $i
      
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