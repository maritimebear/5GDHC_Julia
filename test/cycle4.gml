graph [
    comment "Directed cycle with 4 nodes"
    directed 1
    label "cycle4"
    node [
        id 1
        type "fixed"
        fixed [
            pressure 101325.0
            temperature 298.15
        ]
    ]
    node [
        id 2
        type "junction"
    ]
    node [
        id 3
        type "junction"
    ]
    node [
        id 4
        type "junction"
    ]
    edge [
        source 1
        target 2
        type "prosumer"
        prosumer [
            massflow 1.0
            deltaT 10.0
        ]
    ]
    edge [
        source 2
        target 3
        type "pipe"
        pipe [
            length 1.0
            diameter 1.0
            dx 0.1
        ]
    ]
    edge [
        source 3
        target 4
        type "pipe"
        pipe [
            length 1.0
            diameter 1.0
            dx 0.1
        ]
    ]
    edge [
        source 4
        target 1
        type "pipe"
        pipe [
            length 1.0
            diameter 1.0
            dx 0.1
        ]
    ]
]

