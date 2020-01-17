dist=`awk -F= '/^NAME/{print $2}' /etc/os-release`;
if [ "$dist" == "Ubuntu" ]; then
  echo "ubuntu"
  sudo apt update -y
  DEBIAN_FRONTEND=noninteractive apt-get upgrade -yq
  # sudo apt upgrade -y
  sudo apt install -y build-essential
  sudo apt-get install -y libsqlite3-dev
  sudo apt install -y libgl1-mesa-dev
  sudo apt  install -y awscli 
  sudo apt install -y emacs
  sudo apt install -y unzip
  sudo apt install -y zip
  sudo apt-get install -y zlib1g-dev 
  sudo apt-get install -y scons 
fi
if [ "$dist" == "Amazon Linux AMI" ]; then
    echo "Amazon Linux AMI";
    sudo yum update -y
    sudo yum install awscli emacs unzip zip scons -y
fi



