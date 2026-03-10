output_dir="../randomPatterns"
mkdir -p "$output_dir"

for length in 8 16 32 64; do
    file="${output_dir}/patterns_${length}.txt"
    : > "$file"   # clear file

    for i in $(seq 1 100); do
        tr -dc 'ACGTN' < /dev/urandom | head -c $length >> "$file"

        if [ "$i" -lt 100 ]; then
            echo >> "$file"
        fi
    done
done