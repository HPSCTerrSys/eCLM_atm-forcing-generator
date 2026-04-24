#!/bin/sh
#
# Wrapper to download ERA5 data for a range of months using
# download_ERA5_input.py.
#
# Usage:
#   ./download_ERA5_input_wrapper.sh [start_date=<yyyy-mm>] [end_date=<yyyy-mm>] \
#       [out_dir=<dir>] [request=<custom_request_file>]
#
# Options:
#   start_date   First month to download (inclusive), format yyyy-mm.
#                Default: 2017-07
#   end_date     First month NOT downloaded (exclusive), format yyyy-mm.
#                Default: 2018-08
#   out_dir      Output directory for downloaded files. Default: cdsapidwn
#   request      Path to a custom CDS API request file. When provided, passed
#                to download_ERA5_input.py via --request. If omitted, the
#                default request defined in download_ERA5_input.py is used.
#
# Prerequisites:
#   CDSAPI must be installed and configured with a user-specific API access
#   token before running this script (see README or era5-forcing docs).
#   Must be executed on a login node (internet access required).
#
set -eo pipefail

# Settings
start_date="2017-07" # yyyy-mm
end_date="2018-08"   # yyyy-mm + 1
out_dir="cdsapidwn"
request=""

# Function to parse input
parse_arguments() {
    for arg in "$@"; do
        key="${arg%%=*}"
        value="${arg#*=}"

        case "$key" in
            start_date) start_date="$value" ;;
            end_date) end_date="$value" ;;
            out_dir) out_dir="$value" ;;
            request) request="$value" ;;
            *) echo "Warning: Unknown parameter: $key" ;;
        esac
    done
}

# Call the function to parse the input arguments
# Users needs to make sure for consistent input
parse_arguments "$@"


# create output directory
mkdir -p $out_dir

# loop over months
current_date=$start_date
while [ "$current_date" \< "$end_date" ]; do
    echo "Processing month: $current_date"

    year="${current_date%%-*}"
    month="${current_date#*-}"

    # start download script with data request
    request_opt=""
    [ -n "$request" ] && request_opt="--request $request"
    ./download_ERA5_input.py --year $year --month $month --dirout $out_dir $request_opt

    # Increment the month, arbitrarily setting unimportant day of month to 1
    # POSIX.1-2024 prescribes that months start at zero and years are since 1900
    current_date=$(perl -MPOSIX -e "print strftime( '%Y-%m', 0, 0, 0, 1,
                                    (split(/-/, '$current_date'))[1],
                                    (split(/-/, '$current_date'))[0]-1900);")
done
