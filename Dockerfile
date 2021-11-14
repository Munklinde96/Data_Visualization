FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
SHELL [ "/bin/bash", "-c"]

RUN apt-get update \
    && apt-get install -y git python3-pip nodejs npm zlib1g sudo curl xorg openbox libnss3 libasound2 libatk-adaptor libgtk-3-0
RUN python3 -m pip install --upgrade pip \
    && python3 -m pip install --upgrade Pillow
RUN pip3 install matplotlib seaborn sklearn python-dotenv plotly pandas requests numpy flask flask_cors \
    && pip3 freeze > requirements.txt
RUN curl -sL https://deb.nodesource.com/setup_14.x | sudo -E bash - \
    && apt-get install -y nodejs

WORKDIR /workspaces/Data_Visualization/my-electron-app

RUN npm install --save-dev electron

ARG USERNAME=vscode
# On Linux, replace with your actual UID, GID if not the default 1000
ARG USER_UID=1000
ARG USER_GID=$USER_UID

RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME -s /bin/bash \
    && mkdir -p /home/$USERNAME/.vscode-server /home/$USERNAME/.vscode-server-insiders \
    && chown ${USER_UID}:${USER_GID} /home/$USERNAME/.vscode-server* \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

USER $USERNAME