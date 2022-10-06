#! /bin/sh

SRC=$1
TRG=$2

curl -k -T ${SRC} -u "wnkxtBG28BvvfUq:" -H 'X-Requested-With: XMLHttpRequest' \https://owncloud.gwdg.de/public.php/webdav/${TRG}
