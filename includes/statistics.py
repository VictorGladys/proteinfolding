def mean(num):
    total = 0
    
    for i in num:
        total += i
        
    return total / len(num)

def stddev(num):
    total = 0
    
    for i in num:
        total += (i - mean(num)) ** 2
        
    return (total / (len(num) - 1)) ** 0.5

def skewness(num):
    # Skewness als maat voor moeilijkheid proteine
    sumnum = 0
    sumden = 0
    
    for i in num:
        sumnum += (i - mean(num)) ** 3
        sumden += (i - mean(num)) ** 2
    numer = sumnum
    denom = (len(num) - 1) * stddev(num) ** 3

    return (numer / denom)

if __name__ == '__main__':
    out = []
    for i in range(10):
        inp = input("Number: ")
        out.append(float(inp))
    print("Mean: " + str(mean(out)))
    print("Standard Deviation: " + str(stddev(out)))
    print("Skewness: " + str(skewness(out)))


        
        
