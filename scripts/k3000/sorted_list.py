import sys

def unique(s): # from http://code.activestate.com/recipes/52560-remove-duplicates-from-a-sequence/
    """Return a list of the elements in s, but without duplicates.

    For example, unique([1,2,3,1,2,3]) is some permutation of [1,2,3],
    unique("abcabc") some permutation of ["a", "b", "c"], and
    unique(([1, 2], [2, 3], [1, 2])) some permutation of
    [[2, 3], [1, 2]].

    For best speed, all sequence elements should be hashable.  Then
    unique() will usually work in linear time.

    If not possible, the sequence elements should enjoy a total
    ordering, and if list(s).sort() doesn't raise TypeError it's
    assumed that they do enjoy a total ordering.  Then unique() will
    usually work in O(N*log2(N)) time.

    If that's not possible either, the sequence elements must support
    equality-testing.  Then unique() will usually work in quadratic
    time.
    """

    n = len(s)
    if n == 0:
        return []

    # Try using a dict first, as that's the fastest and will usually
    # work.  If it doesn't work, it will usually fail quickly, so it
    # usually doesn't cost much to *try* it.  It requires that all the
    # sequence elements be hashable, and support equality comparison.
    u = {}
    try:
        for x in s:
            u[x] = 1
    except TypeError:
        del u  # move on to the next method
    else:
        return u.keys()

    # We can't hash all the elements.  Second fastest is to sort,
    # which brings the equal elements together; then duplicates are
    # easy to weed out in a single pass.
    # NOTE:  Python's list.sort() was designed to be efficient in the
    # presence of many duplicate elements.  This isn't true of all
    # sort functions in all languages or libraries, so this approach
    # is more effective in Python than it may be elsewhere.
    try:
        t = list(s)
        t.sort()
    except TypeError:
        del t  # move on to the next method
    else:
        assert n > 0
        last = t[0]
        lasti = i = 1
        while i < n:
            if t[i] != last:
                t[lasti] = last = t[i]
                lasti += 1
            i += 1
        return t[:lasti]

    # Brute force is all that's left.
    u = []
    for x in s:
        if x not in u:
            u.append(x)
    return u

def compare (tuple1,tuple2):
    '''
    if tuple1 starts with tuple2: return 0
    if tuple1<tuple2: return -1
    if tuple1>tuple2: return 1
    '''
    if len(tuple1)>0 and type(tuple1[-1]) is str: tuple1=tuple1[:-1] # TODO: optimization feasible in calling functions for avoiding those systematic tests
    if len(tuple2)>0 and type(tuple2[-1]) is str: tuple2=tuple2[:-1] # TODO: optimization feasible in calling functions for avoiding those systematic tests
    tmp_tuple1=tuple1[0:len(tuple2)]
    if tmp_tuple1 < tuple2: return    -1
    if tmp_tuple1 > tuple2: return     1
    return                             0

class sorted_list(object):
    """Class sorted list
    Contains a sorted set of list of integers (positives or negatives, not equal to 0)
    Divided into buckets. Each first value is a bucket.
    """

    def __init__(self):
        self.main_dict={}
        self.size=0
        

    def add(self, mylist):
        """add a new list"""
        self.size+=1
        if mylist[0] not in self.main_dict: 
            self.main_dict[mylist[0]]=[]
        self.main_dict[mylist[0]]+=[mylist[1:]]
        
    
    def sorted_add(self,mylist):
        """add a new list"""
        self.add(mylist)
        
        # we added the element at the last position of the list. We will bring it back to the good position:
        # example the list is: 4,N,7,13,N,N,5 (we have to put back 5 in its good location). Note that in practice we 
        # pos_current = 6, previous = 5, 4, 3.  current_list[3]<current_list[6]: yes, so we swap, and obtain 4,N,7,5,N,N,13
        # pos_current = 3, previous = 2         current_list[2]<current_list[3]: yes, so we swap, and obtain 4,N,5,7,N,N,13
        # pos_current = 2, previous = 1, 0      current_list[0]<current_list[2]; no, this is finished. 
        current_list = self.main_dict[mylist[0]]
        pos_current = len(current_list)-1
        while True:
            previous = pos_current-1
            while current_list[previous]==None and previous>=0:
                previous-=1
                
            if previous==-1: break                                          # end of the list, nothing to do.
            
            if current_list[previous] <= current_list[pos_current] : break  # bubble sort of the last element finished. Note that the "equal case" should not happen in our situation. 
            else :                                                          # we swap the two values
                current_list[previous], current_list[pos_current] = current_list[pos_current], current_list[previous]
                pos_current     =   previous;


        

    def sort(self):
        """sort the whole structure - feasible only before removing elements"""
        for key, value in self.main_dict.items():
            value.sort()
            
    def traverse(self):
        for key, value in self.main_dict.items():
            # mylists = value
            for mylist in value: 
                if mylist != None:
                    yield [key]+mylist
                    
    def index_nodes(self):
        ''' For each element in the structure, we add its id as a last value, stored as i_14, for instance for node 14.'''
        index_id=0
        for key, value in self.main_dict.items():
            for mylist in value: 
                if mylist != None:
                    mylist+=['i_'+str(index_id)]
                    index_id+=1
                    
    
    def remove(self,mylist):
        '''remove an element from the structure'''
        ''' if the element is not in a structure a warning is raised'''
        try:
            tormindex=self.main_dict[mylist[0]].index(mylist[1:])
            del self.main_dict[mylist[0]][tormindex]
        except :
            sys.stderr.write("\n        WARNING: Tried to remove "+str(mylist)+" absent from the list.\n")
            return 0
        self.size-=1
        return 1

    def get_lists_starting_with_given_prefix(self, prefix):
        ''' given a prefix of a list, return all lists in the set starting with this prefix
        Dichotomy. log(|self|/size(alphabet)) comparisons
        O(|prefix|*log(|self|/size(alphabet)))
        '''
        first=prefix[0]
        if first not in self.main_dict: return []
        current_list = self.main_dict[first]
        prefix = prefix[1:]
        start=0
        n=len(current_list)-1
        stop=n
    
   
        # print ("search", prefix, "in ",current_list)
        while start <=stop: 
            # print (start,stop)
            middle = start+int((stop-start)/2)
            if middle>n: break
            ## Search for a non empty middle value
            while middle<=stop:
                tuple1 = current_list[middle]
                if tuple1 == None:   middle+=1  # Empty value
                else:                break      # Found a non empty value
            if middle == stop+1:                 # for instance looking for [1] in [1, None, None] we didn't find any middle with a avlue, thus we come back 
                middle = start+int((stop-start)/2)-1
                while middle>=start:
                    tuple1 = current_list[middle]
                    if tuple1 == None:   middle-=1  # Empty value
                    else:                break      # Found a non empty value
            if middle == start-1: return []     # nothing found (empty sub-list)
                
            ## a non empty middle value was found
            cmp_val = compare(tuple1,prefix)
            if cmp_val == -1:   # prefix may start in the remaining part of the current_list array
                start = middle+1
                continue
            if cmp_val ==  1:   # prefix may start in the starting part of the current_list array
                stop = middle-1
                continue
            # if cmp_val == 0:    # we found a tuple starting with the prefix. We need to check other tuples starting with prefix before and after in the array.
            res=[[first]+tuple1]
            # print ("res before all = ",res)
            i=middle-1
            while i>-1:
                if current_list[i]!=None:
                    if compare(current_list[i],prefix)==0: res.append([first]+current_list[i])
                    else: break
                i-=1
                
            i=middle+1
            while i<=n: 
                if current_list[i]!=None:
                    if compare(current_list[i],prefix)==0: res.append([first]+current_list[i])
                    else: break
                i+=1
                    
            # print ("res after all = ",res)
            return res
        return []
        
    def contains (self,mylist):
        ''' if sr is in SR: return True, else return False.
        faster (O(log(|SR|)) that the 'in' function from non ordered array (O(|SR|))
        '''
        X=self.get_lists_starting_with_given_prefix(mylist) 
        if X==mylist: return True   # if X is composed of a unique list
        if mylist in X: return True # O(|X|) which can be considered as constant
        return False
        
    def __len__(self):
        return self.size
        
        
    def contains(self,query):
        """only for debug purpose, no need to be optimized"""
        for key, value in self.main_dict.items():
            if key==query: 
                print("key", key, "query",query) 
                return True
            for followup in value:
                for unitig_id in followup:
                    if unitig_id==query: 
                        return True
        return False
        
    def unique(self):
        for key, value in self.main_dict.items():
            size_before=len(value)
            value=unique(value)
            size_after=len(value)
            self.main_dict[key]=value
            self.size-=(size_before-size_after)
            
    def __str__(self):
        str_SR=""
        for key, value in self.main_dict.items():
            for followup in value:
                str_SR+=str(key)
                for unitig_id in followup:
                    str_SR+=","+str(unitig_id)
                str_SR+="\n"
        return str_SR

# SR = sorted_list()
# SR.add([1,2,3])
# SR.add([1,2,3])
# SR.remove([1,2])
# # SR.remove([1])
# # SR.remove([1,2,4,5])
# SR.add([1,2,3,5])
# # SR.add([5])
# # SR.add([5,6])
# print(SR)
# SR.unique()
# print(SR)