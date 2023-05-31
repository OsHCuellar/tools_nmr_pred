#!/bin/bash
OLD_INSTANCE=`docker ps -q -f name=tools-example-instance`
if [ "$OLD_INSTANCE" != "" ]
then
    docker kill $OLD_INSTANCE
fi
docker run -d -p 8091:80 --rm --name=tools-example-instance tools-example && echo "You can connect to http://localhost:8091"