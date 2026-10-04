using Dates

function current_unix_time()
    return Dates.datetime2unix(Dates.now())
end

function elapsed_time(start_time)
    return current_unix_time() - start_time
end

# Example usage
print(Dates.now())
start_time = current_unix_time()
print("tiempo actual $start_time\n")

# Simulate waiting for some time
sleep(2)  # Sleep for 2 seconds
delta = elapsed_time(start_time)
print("Elapsed time:  $delta seconds")
