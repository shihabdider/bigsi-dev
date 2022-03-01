function main() {
    const bin = 11
    console.log(bin.toString(2))
    const query = 12
    console.log(query.toString(2))
    const intersection = bin & query
    console.log(intersection, intersection.toString(2))
}

main()
