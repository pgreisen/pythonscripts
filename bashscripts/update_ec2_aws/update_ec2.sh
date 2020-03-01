dist=`awk -F= '/^NAME/{print $2}' /etc/os-release`;
if [ "$dist" == "\"Ubuntu\"" ]; then
  echo "ubuntu"
  sudo apt update -y;
  sudo DEBIAN_FRONTEND=noninteractive apt-get upgrade -yq;
  sudo apt install -y build-essential;
  sudo apt-get install -y libsqlite3-dev;
  sudo apt install -y libgl1-mesa-dev;
  sudo apt  install -y awscli python3-pip emacs unzip zip scons
  sudo apt-get install -y zlib1g-dev 
fi
if [ "$dist" == "\"Amazon Linux\"" ]; then
    echo "Amazon Linux AMI";
    sudo yum update -y
    sudo yum install awscli emacs unzip zip scons python3-pip -y
fi



