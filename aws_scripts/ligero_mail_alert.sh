#!/bin/bash

while IFS=, read -r email
    do
       echo $email
    echo "Report attached." | mail -s "$ALERT_MAIL" "$email" -A AWS/usage.csv
    done < AWS/mail_alerts
