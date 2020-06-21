dist=`awk -F= '/^NAME/{print $2}' /etc/os-release`;
if [ "$dist" == "\"Ubuntu\"" ]; then
  echo "ubuntu"
  sudo apt-get update -y;
  sudo DEBIAN_FRONTEND=noninteractive apt-get upgrade -yq;
  sudo apt-get install -y build-essential;
  sudo apt-get install -y libsqlite3-dev;
  sudo apt-get install -y libgl1-mesa-dev;
  sudo apt-get install -y awscli python3-pip emacs unzip zip 
  sudo apt-get install -y zlib1g-dev 
  sudo apt-get install -y scons  
fi
if [ "$dist" == "\"Amazon Linux\"" ]; then
    echo "Amazon Linux AMI";
    sudo yum update -y
    sudo yum install awscli emacs unzip zip scons python3-pip -y
fi
