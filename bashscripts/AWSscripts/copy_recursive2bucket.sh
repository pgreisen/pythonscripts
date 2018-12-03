aws s3 cp s3://data/ . --recursive --exclude "*" --include "2016-08*"`

# Region can be necessary to specify
# aws s3 cp s3://p450/PSSM . --recursive --include "*" --region eu-central-1