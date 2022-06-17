function test_array() {
    local msg="$1"   # Save first argument in a variable
    shift            # Shift all arguments to the left (original $1 gets lost)
    local arr=("$@")
    for i in ${arr[@]}; do echo "$i"; done
}

arr=(1 2 3)
test_array hello ${arr[@]}
