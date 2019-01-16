echo "docker run --rm -it -v $PWD:/pdx-hts -p 8888:8888 -e GRANT_SUDO=yes edjuaro/pdx-hts:3.0 /pdx-hts/launch_jupyter.sh"
docker run --rm -it -v $PWD:/pdx-hts -p 8888:8888 -e GRANT_SUDO=yes edjuaro/pdx-hts:3.0 /pdx-hts/launch_jupyter.sh
