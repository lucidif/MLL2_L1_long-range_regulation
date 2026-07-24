#!/usr/bin/env bash
#set -euo pipefail

HOST_RAW="$(hostname 2>/dev/null || cat /etc/hostname)"
HOST="${HOST_RAW%%.*}"
HOST="${HOST,,}"

# Default
DOCKER_ARGS=(
  "-u" "$(id -u):$(id -g)"
  "-v" "${PWD}:${PWD}"
  "-w" "${PWD}"
)

# Override per host
case "$HOST" in
  ziggystardust)
    DOCKER_ARGS=(
      "-u" "$(id -u):$(id -g)"
      "-v" "${PWD}:${PWD}"
      "-v" "/mnt/datawk1:/mnt/datawk1"
      "-v" "/media/lucio/easystore:/media/lucio/easystore"
      "-w" "${PWD}"
    )
    ;;
  travelmatto)
    DOCKER_ARGS=(
      "-u" "$(id -u):$(id -g)"
      "-v" "${PWD}:${PWD}"
      "-v" "/mnt/bowie/cache:/mnt/bowie/cache"
      "-w" "${PWD}"
    )
    ;;
esac

export HOST

#how to use it
#source "git/Lara_MLL2/bin/docker_env.sh"
#sudo docker run --rm "${DOCKER_ARGS[@]}" -it quay.io/biocontainers/pairix:0.3.7--py36h30a8e3e_3