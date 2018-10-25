'''
aws configure
AWS Access Key ID [None]: XXXXXXXXXXXXXXXXXXX
AWS Secret Access Key [None]: XXXXXXXXXXXXXXXXX
Default region name [None]:  enter
Default output format [None]:  enter
'''
file2transfer=abc.tgz
path2transferto=mypath/mypath
aws s3 cp $file2transfer s3://$path2transferto