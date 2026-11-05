# Commands written by ADIT. ADIT does not run them; paste them into a terminal.
# Replace the <...> parts with your own values (set host and remote_dir in the profile to fill them in).

# Work out where this file is, so the commands keep working if the folder is moved
HERE="$(cd "$(dirname "$0")" && pwd)"

# 1. Send this directory to the cluster
rsync -av "$HERE" <user>@<cluster host name>:<work directory on the cluster>/

# 2. Log in to the cluster and submit
ssh <user>@<cluster host name>
cd <work directory on the cluster>/vasp_cd_si_generated && qsub submit.sh
qstat -u $USER

# 3. Bring the results back when it has finished
rsync -av <user>@<cluster host name>:<work directory on the cluster>/vasp_cd_si_generated "$(dirname "$HERE")"/
