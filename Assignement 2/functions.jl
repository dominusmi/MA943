# Recurrent function which trasverses a linked list and prints the data of each node
# @param:  list::Nullable{LList}
# @return: none
function list_trasverser(list::Nullable{LList})
    if isnull( list ) 
        # Check whether at leaf node
        return;
    else
        # If not, show data 
        node = get(list)
        @show node.data
        list_trasverser(node.next)
    end
end

# 'Linearly' searches a linked list for the node with a key value 'k'
# @param 
# list: linked list to search
# k: key to search for 
# @return 
# ::Pair: Pairnode with key 'k' if found, othw 'false'
function search_list(list::Nullable{LList}, k::Int)
    if isnull( get(list).next )
        # Check whether at leaf node
        return false;
    else
        # If not, check current node key, and otherwise pass on
        node = get(list)
        if( node.data.key == k )
            return node.data
        else
            return search_list(node.next, k)
        end
    end  
end

# Generates a linked list of random number of size 'size', with increasing unique keys in 1:'size'
# @param 
# size: size of the linked list
# @return 
# list::Nullable{LList} 
function generate_llist(size::Int)
    arr = Array(KVPair,0)
    for i in 1:size
        tmp = KVPair(i,rand())
        push!(arr, tmp)
    end

    return buildLList(arr)
end

# Generates a linked list containing the cumulative sum
# @param: 
# size: sisze of list to create
# rng: random generator
# @return
# Nullable{LList}: The generated list
# Array{floats}: the array containing the cumulative sum at every node
function generate_cumulative_llist(size::Int, rng)
    # Set number of intervals
    n=size

    # Generate a sample of n uniform random numbers in the interval [0,1]
    X = rand(rng, n);

    # Generate array of key-values
    pairs = Array(KVPair,0)
    cumsum = zeros(n+1)
    for i in 1:n
        cumsum[i+1] = cumsum[i] + X[i]
        push!(pairs, KVPair(i, cumsum[i+1]))

    end
    
    return buildLList(pairs), cumsum
end

# 'Linearly' searches linked list for the node representing the interval containing 'x'
# @params
# list: list to search
# x: float whose interval is being searched
# @return
# Pair: pair of key-value which contains x
function intervalmembership_brute(list::Nullable{LList}, x::Float64)
    if isnull(list)
        return -1
    else
        node = get(list)
        if node.data.value > x
            return node.data
        else
            return intervalmembership_brute(node.next,x)
        end
    end
end

# Searches the Fenwick tree for the node representing the interval containing 'x'
# @param
# FT: Tree to search
# x: value to search
# @return
# Pair: pair object containing the value x
function intervalmembership(FT::Nullable{FTree}, x::Float64)
    if isnull(FT)
        # Really shouldn't happen ...
        return -1
    else
        # Normal case 
        node = get(FT)
        if isnull( node.right ) && isnull( node.left )
            # If at leaf, return current
            return node.data
        else
            # Else check which tree to return
            y = get(node.left).data.value
            if x <= y
                return intervalmembership(node.left, x)
            else
                return intervalmembership(node.right, x-y)
            end
        end
    end     
end

# Generates intervals and cumulative sums -- only used for testing
function generate_intervals(size, rng)
    
    n=size
    X = rand(rng, n)

    # Now calculate the array of partial sums
    values = Array{KVPair}(n)
    cumsum = zeros(n+1)

    for i in 1:n
        values[i] = KVPair(i,X[i])
        cumsum[i+1] = cumsum[i]+X[i]
    end
    
    return values,cumsum,rng
end

