echo "docker run --rm -v $PWD:/pdx-hts -p 8888:8888 -e GRANT_SUDO=yes edjuaro/pdx-hts:2.0 /pdx-hts/launch_jupyter.sh"
docker run --rm -v $PWD:/pdx-hts -p 8888:8888 -e GRANT_SUDO=yes edjuaro/pdx-hts:2.0 /pdx-hts/launch_jupyter.sh
