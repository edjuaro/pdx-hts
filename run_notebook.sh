echo "docker run --rm -v $PWD/pdx-hts:/pdx-hts -p 8888:8888 -e GRANT_SUDO=yes edjuaro/pdx-hts:1.0 /pdx-hts/launch_jupyter.sh"
# To run a notebook
# docker run --rm -v $PWD/mb-pdx-hts:/mb-pdx-hts -p 8888:8888 -e GRANT_SUDO=yes edjuaro/molecular-tumor-board:6.0 /mb-pdx-hts/launch_jupyter.sh
docker run --rm -v $PWD/pdx-hts:/pdx-hts -p 8888:8888 -e GRANT_SUDO=yes edjuaro/pdx-hts:1.0 /pdx-hts/launch_jupyter.sh
