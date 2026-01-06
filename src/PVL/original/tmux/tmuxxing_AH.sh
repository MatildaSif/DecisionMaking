#!/usr/bin/env bash

##### HUGE KUDOS TO Balázs FOR SETTING THIS UP!!!! #####

#install packages
source setup.sh

SESSION_NAME="TMUXXING"

#commands=("tail -f log1" "tail -f log2" "tail -f log3")
command="tmux send-keys -t $i 'Rscript Tmux_ORL_rec.R '"
#tmux_command+=" ${commands[0]}"
#for cmd in "${commands[@]:1}"; do
#    tmux_command+="\; split-window $cmd"
#done
#tmux "$SESSION_NAME \; attach"

# Check if the session already exists
if tmux has-session -t $SESSION_NAME 2>/dev/null; then
    echo "Session $SESSION_NAME already exists. Attaching to it."
    tmux attach-session -t $SESSION_NAME
else
    # Create a new session and panes 
    # (there must be a more automated way with a loop)
    #tmux new-session -d -s $SESSION_NAME
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
    
    SEED=1000
    #tmux_command+=" tmux send-keys -t 0 'Rscript Tmux_ORL_rec.R ' $SEED C-m"
    tmux send-keys -t 0 'Rscript Tmux_ORL_rec.R ' $SEED C-m
    
    #loop for seeds and send them to the script
    for i in {1..20}
    do
      #define seed
      SEED=$((i + 1001))
      
      # build commands
      #tmux_command+="\; split-window tmux send-keys -t $i 'Rscript Tmux_ORL_rec.R ' $SEED C-m"
      tmux split-window send-keys -t $i 'Rscript Tmux_ORL_rec.R ' $SEED C-m
      #echo $tmux_command
      
      # Send a command to the pane, make R read SEED
      #tmux send-keys -t $i 'Rscript Tmux_ORL_rec.R ' $SEED C-m
    
    done
    tmux "$SESSION_NAME \; attach"
 # Attach to the created session
  #tmux attach-session -t $SESSION_NAME
fi