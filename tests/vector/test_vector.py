import pytest

from protein_folding.vector import Vec3D


def test_create_vec():
    _ = Vec3D(0, 0, 0)


def test_abs():
    v = Vec3D(1, -4, 5)
    v_abs = abs(v)

    assert all((v_abs.x == 1, v_abs.y == 4, v_abs.z == 5))


def test_add():
    v1 = Vec3D(1, 4, 2)
    v2 = Vec3D(4, 2, 1)
    v = v1 + v2

    assert all((v.x == 5, v.y == 6, v.z == 3))


def test_iadd():
    v1 = Vec3D(1, 4, 2)
    v2 = Vec3D(4, 2, 1)
    v1 += v2

    assert all((v1.x == 5, v1.y == 6, v1.z == 3))


def test_iter():
    v = Vec3D(1, 4, 2)
    x, y, z = v

    assert all((x == 1, y == 4, z == 2))


def test_sub():
    v1 = Vec3D(1, 4, 2)
    v2 = Vec3D(4, 2, 1)
    v = v1 - v2

    assert all((v.x == -3, v.y == 2, v.z == 1))


def test_abs_diff():
    v1 = Vec3D(1, 4, 2)
    v2 = Vec3D(4, 2, 1)
    v = Vec3D.abs_diff(v1, v2)

    assert all((v.x == 3, v.y == 2, v.z == 1))


def test_area():
    v = Vec3D(1, 4, 2)
    v_area = v.area()

    assert v_area == 4


def test_eq_components():
    v = Vec3D(1, 4, 2)
    eq = v.eq_components(2)

    assert eq == (False, False, True)


def sum_components():
    v = Vec3D(1, 4, 2)
    s = v.sum_components()

    assert s == 7


def test_len_sq():
    v = Vec3D(1, 4, 2)
    l = v.len_sq()

    assert l == 21


def test_volume():
    v = Vec3D(1, 4, 2)
    volume = v.volume()

    assert volume == 8
