#Backs up the entire spring_2022_study folder to keoki, and then Box.

set -e #exit the script if there are errors; this guarantees that both syncs were successful if the script runs through

rsync -av --delete /Users/dp/projects/ARL/ /Volumes/experiments2/Daniel/ARL #syncs to keoki, deleting any files on keoki that aren't on your machine (you can access old files on Box)
rclone copy /Users/dp/projects/ARL Duke_box:/ARL -v #copies to Box
