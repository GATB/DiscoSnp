/********************************************************************************/
/* Using Tool class for quick tool development. */
/********************************************************************************/
// We define some constant strings for names of command line parameters
static const char* STR_RANGE_MIN = "-min";
static const char* STR_RANGE_MAX = "-max";
/********************************************************************************/
// We define our own tool as a class that inherits from the Tool class
// By doing this, we will get:
// - a command line parser with default options; we will add specific options
// - a mean to retrieve command line options values in our code
// - a mean to get configurable progression bar for our ongoing job
// - a dispatcher object for easy multi-threading
// - a mean to get time execution information
// - a mean to gather any piece of information for output
// - a mean to run the tool in the 'main' function

    // The Tool constructor allows to give a name to our tool.
    // This name appears when one gets help in the command line or in the final output
    Kissreads2 () : Tool ("Kissreads2")
    {
        // By default, the Tool class provides 3 default command line arguments:
        // "-nb-cores X" : number of cores used by the dispatcher object (0 means all cores)
        // "-verbose X" : verbosity level (0: none, 1 and 2: bargraphs and final output)
        // "-help" : just dump help about the tool
        // Now, we add two custom command line arguments with the parser we got from Tool
        // the "-min" argument is not mandatory, with default value 1
        // the "-max" argument is mandatory
        getParser()->push_front (new OptionOneParam (STR_RANGE_MIN, "lower range bound", false, "1"));
        getParser()->push_front (new OptionOneParam (STR_RANGE_MAX, "upper range bound", true));
        // Hence our tool can get 5 arguments as input (3 default, 2 custom)
        // One can get them by using "./Kissreads2 -help"
        // Note: if the mandatory argument is not provided, an error is dumped in console
    }
    
    
    // The 'execute' method must be implemented as we are a subclass of Tool
    // This method does the actual job our tool is supposed to do (here some dummy computation)
    void execute ()
    {
        // Our job is to compute some dummy information from integers in an range.
        // This range is configured with our two command line arguments.
        // Here, we see how to retrieve the arguments values through the 'getInput' method
        Range<u_int64_t> range (
                                getInput()->getInt(STR_RANGE_MIN),
                                getInput()->getInt(STR_RANGE_MAX)
                                );
        // We create an iterator over our integer range.
        // Note how we use the Tool::createIterator method. According to the value of the "-verbose" argument,
        // this method will add some progression bar if needed.
        Iterator<u_int64_t>* iter = createIterator<u_int64_t> (range, "iterate range");
        LOCAL (iter);
        // We will do some dummy computation
        u_int64_t totalSum = 0;
        // We want to get execution time. We use the Tool::getTimeInfo() method for this.
        getTimeInfo().start ("computation");
        // We iterate the range through the Dispatcher we got from our Tool parent class.
        // The dispatcher is configured with the number of cores provided by the "-nb-cores" command line argument.
        IDispatcher::Status status = getDispatcher()->iterate (iter, [&] (const u_int64_t& i)
                                                               {
                                                                   cout << "i="<<i<<endl;
                                                                   // We do some dummy computation.
                                                                   u_int64_t sum = 0;
                                                                   for (u_int64_t j=0; j<i; j++) { sum += j; }
                                                                   __sync_fetch_and_add (&totalSum, sum);
                                                               });
        getTimeInfo().stop ("computation");
        // We gather some statistics. Note how we use the getInfo() method from the parent class Tool
        // If the verbosity is not 0, all this information will be dumped in the console at the end
        // of the tool execution
        getInfo()->add (1, "output");
        getInfo()->add (2, "sum", "%ld", totalSum);
        getInfo()->add (1, getTimeInfo().getProperties("time"));
    }
