#!/bin/bash
# Open a notebook, first run on geir/yngve, then on lanermac

# Where do you call Jupyter from on the server?
# JUPYTER_LOC=jupyter
REMOTE_JUPYTER_LOC=~/.local/bin/jupyter

# Set the ports
read -p 'Remote port [8897] ' REMOTE_PORT 
REMOTE_PORT=${REMOTE_PORT:-8897}
read -p 'Local port [8889] ' LOCAL_PORT
LOCAL_PORT=${LOCAL_PORT:-8889}

# Set the machine names
GEIR_MACHINE=geir.astro.utoronto.ca
YNGVE_MACHINE=yngve.astro.utoronto.ca
LOCAL_MACHINE=lanermac

# If we're on the server then open the notebook
if [ $HOSTNAME == geir ] || [ $HOSTNAME == geir.astro.utoronto.ca ] || [ $HOSTNAME == yngve ] || [ $HOSTNAME == yngve.astro.utoronto.ca ]; then
  read -p 'jupyter type [*lab/notebook] ' JUPYTER_TYPE
  JUPYTER_TYPE=${JUPYTER_TYPE:-lab} # Default is lab
  read -p 'launch directory [../] ' RUN_DIR
  RUN_DIR=${RUN_DIR:-../} # Default is root project directory
  echo 'Opening a jupyter '$JUPYTER_TYPE' session in '$RUN_DIR
  $REMOTE_JUPYTER_LOC $JUPYTER_TYPE --no-browser --port=$REMOTE_PORT $RUN_DIR
fi

# Open the SSH tunnel
if [ $HOSTNAME == $LOCAL_MACHINE ] || [ $HOSTNAME == $LOCAL_MACHINE.local ] || [ $HOSTNAME == $LOCAL_MACHINE.lan ] ; then
  read -p 'server name: [*geir/yngve] ' SERVER
  SERVER=${SERVER:-geir} # Default is geir
  echo 'Tunnel open to '$SERVER' through local port '$LOCAL_PORT
  ssh -N -L localhost:$LOCAL_PORT:localhost:$REMOTE_PORT $SERVER
fi
