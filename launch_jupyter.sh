#!/bin/sh
echo "->We are now running the analysis."
cd /pdx-hts/Notebooks

echo "You may have dropped this from your clipboard:"
echo "\t$pbpaste"
echo "USE THIS URL:"
echo "\t\thttp://127.0.0.1:8888/"
echo "http://127.0.0.1:8888/" | pbcopy
jupyter-notebook --no-browser --port 8888 --ip=0.0.0.0 --NotebookApp.password='' --NotebookApp.token='' --NotebookApp.password_required=False --allow-root
