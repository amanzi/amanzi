#!/usr/bin/env bash

mkdir ~/.jupyter 
jupyter_cfg=~/.jupyter/jupyter_notebook_config.py
echo "c.JupyterApp.config_file = ''" >> $jupyter_cfg 
echo "c.NotebookApp.allow_root = True" >> $jupyter_cfg 
echo "c.NotebookApp.allow_remote_access = True" >> $jupyter_cfg
echo "c.NotebookApp.ip = '*'" >> $jupyter_cfg
echo "c.NotebookApp.terminado_settings = { \"shell_command\": [\"/usr/bin/bash\"] }" >> $jupyter_cfg
echo "c.NotebookApp.token = u''" >> $jupyter_cfg
