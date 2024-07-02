
n to get the modification time in hours since the epoch
get_mod_time_in_hours() {
  local file=$1
  if [ -e "$file" ]; then
    # Get the modification time in seconds since the epoch
    mod_time=$(stat -c %Y "$file")
    # Convert seconds to hours
    echo "scale=2; $mod_time / 3600" | bc -l
  else
    echo "File named '$file' does not exist."
    exit 1
  fi
}

# Function to count lines containing "ip:" in a file
count_ip_lines() {
  local file=$1
  local ip_count=$(grep -c 'ip:' "$file")
  echo "$ip_count"
}

# Read number of tasks and CPUs per task from run_slurm1.sh
ntasks=$(grep -oP '(?<=^#SBATCH --ntasks=)\d+' run_slurm1.sh)
cpus_per_task=$(grep -oP '(?<=^#SBATCH --cpus-per-task=)\d+' run_slurm1.sh)

# Debugging: Print ntasks and cpus_per_task
echo "ntasks: $ntasks"
echo "cpus_per_task: $cpus_per_task"

if [ -z "$ntasks" ] || [ -z "$cpus_per_task" ]; then
  echo "Failed to read ntasks or cpus-per-task from run_slurm1.sh"
  exit 1
fi

# Calculate nthnum
nthnum=$((ntasks * cpus_per_task))

# Debugging: Print nthnum
echo "nthnum: $nthnum"

# Get the modification time in hours for pes and stdout.log files
pes_hours=$(get_mod_time_in_hours "pes")
stdout_log_hours=$(get_mod_time_in_hours "stdout.log")

# Debugging: Print modification times
echo "pes_hours: $pes_hours"
echo "stdout_log_hours: $stdout_log_hours"

# Calculate the time difference in hours
time_diff=$(echo "scale=2; $stdout_log_hours - $pes_hours" | bc -l)

# Debugging: Print time difference
echo "time_diff: $time_diff"

# Count lines containing "ip:" in stdout.log
ip_count=$(count_ip_lines "stdout.log")

# Debugging: Print ip_count
echo "ip_count: $ip_count"

# Extract nenrg and nkline values from stdout.log
nenrg=$(awk -F ': *' '/nenrg:/ { print $2; found=1; exit } END { if (!found) print "0" }' stdout.log)
nkline=$(awk -F ': *' '/nkline:/ { print $2; found=1; exit } END { if (!found) print "0" }' stdout.log)

# Debugging: Print nenrg and nkline
echo "nenrg: $nenrg"
echo "nkline: $nkline"

# Calculate nppt as nenrg * nkline
nppt=$(echo "scale=2; $nenrg * $nkline" | bc -l)

# Debugging: Print nppt
echo "nppt: $nppt"

# Calculate the required value
if [ "$nthnum" -eq 1 ]; then
  echo "nthnum cannot be 1 to avoid division by zero."
  exit 1
fi

# Calculate result using the formula: ceil(nppt / nthnum) * (time_diff / (ip_count / nthnum - 1))
ceil_nppt=$(echo "($nppt + $nthnum - 1) / $nthnum" | bc)
result=$(echo "scale=2; $ceil_nppt * ($time_diff / ($ip_count / $nthnum - 1))" | bc -l)

# Output the result
echo "Result: $result"

