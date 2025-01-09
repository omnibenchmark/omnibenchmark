onsuccess:
    shell('find out -name "*performance.txt" | sort | xargs head > performances.txt')
