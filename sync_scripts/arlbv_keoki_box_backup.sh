set -e #exit the script if there are errors; this guarantees that both syncs were successful if the script runs through

#rsync -av --delete /Users/dp/projects/ARL_bv/ /Volumes/experiments2/Daniel/ARL_bv #syncs to keoki, deleting any files on keoki that aren't on your machine (you can access old files on Box)
rclone copy /Users/dp/projects/ARL_bv Duke_box:/PROJECT\ 3373\:\ Projects/ARL_bv -v #copies to Box
